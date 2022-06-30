/* Project specific Javascript goes here. */
function populateResponse(response){
  document.getElementById('captcha-response').value = response;
}
function showLoader(){
  $("#loading-overlay").show();
}
function hideLoader(){
  $("#loading-overlay").hide();
}

function showHideLoader(time) {
  showLoader();
  await delay(time);
  hideLoader();
}
function submitForm(form){
  form.submit();
}
function redirectTo(url){
  window.location=url;
}
function showLoaderOnClick(url) {
    showLoader();
    window.location=url;
  }

function showLoaderOnClickNoURL() {
    showLoader();
  }

function submitRecaptcha(form, button){
  button.on('click', function(e){
    var response = $("#captcha-response");
    if (response == ""){
      e.preventDefault();
      alert("Cannot validate ReCAPTCHA. Please try again.");
    } else {
      showLoader();
      form.submit();
    }
  });
}
